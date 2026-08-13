#!/usr/bin/env python3
"""Independently reconstruct Hardy telemetry chains and compare three proof routes."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
FIXTURE_ROOT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10"
)
FIXTURE_RESULT = FIXTURE_ROOT / "fixture_result.json"
TELEMETRY = FIXTURE_ROOT / "enabled/chain_telemetry.jsonl"
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_chain_telemetry_scout.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_c1_hardy_chain_telemetry_scout.md"
)
SCRIPT = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_hardy_chain_telemetry_scout.py"
PRECISION_LADDER = (50, 80, 120)
mp.mp.dps = max(PRECISION_LADDER) + 40


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        require(line.strip(), f"blank telemetry record at line {line_number}")
        records.append(json.loads(line))
    require(records, "empty telemetry file")
    return records


def mpc(value: list[str]) -> mp.mpc:
    require(isinstance(value, list) and len(value) == 2, "invalid complex record")
    result = mp.mpc(mp.mpf(value[0]), mp.mpf(value[1]))
    require(mp.isfinite(result.real) and mp.isfinite(result.imag), "non-finite complex record")
    return result


def decimal(value: mp.mpf | mp.mpc, digits: int = 40) -> str | list[str]:
    if isinstance(value, mp.mpc):
        return [mp.nstr(value.real, digits), mp.nstr(value.imag, digits)]
    return mp.nstr(value, digits)


def norm(value: mp.mpc) -> mp.mpf:
    return abs(value)


def independent_direct_sum(level: dict[str, Any], tpm: mp.mpf, dps: int) -> mp.mpc:
    length = int(level["length"])
    coefficients = [mp.mpf(value) for value in level["coefficients"]]
    require(length >= 0, "negative direct level length")
    require(len(coefficients) == int(level["degree"]), "coefficient/degree mismatch")

    with mp.workdps(dps):
        total = mp.mpc(0)
        for index in range(length + 1):
            y = mp.mpf(index)
            phase = mp.mpf(0)
            for coefficient in reversed(coefficients):
                phase = (phase + coefficient) * y
            total += mp.exp(mp.j * tpm * phase)
        return +total


def adapt(level: dict[str, Any], raw_value: mp.mpc) -> mp.mpc:
    value = mp.conj(raw_value) if level["conjugate"] else raw_value
    if level["subtract_one"]:
        value -= 1
    return value


def apply_step(step: dict[str, Any], child: mp.mpc) -> mp.mpc:
    value = mpc(step["multiplier"]) * child + mpc(step["qq"])
    if step["conjugate"]:
        value = mp.conj(value)
    if step["subtract_one"]:
        value -= 1
    return value


def apply_homogeneous(step: dict[str, Any], error: mp.mpc) -> mp.mpc:
    value = mpc(step["multiplier"]) * error
    return mp.conj(value) if step["conjugate"] else value


def quantile(values: list[mp.mpf], fraction: float) -> mp.mpf:
    require(values, "empty quantile input")
    ordered = sorted(values)
    index = round((len(ordered) - 1) * fraction)
    return ordered[index]


def finite_ratio(numerator: mp.mpf, denominator: mp.mpf) -> mp.mpf | None:
    if denominator == 0:
        return None
    return numerator / denominator


def ratio_text(value: mp.mpf | None) -> str | None:
    return None if value is None else str(decimal(value, 30))


def reduction_text(numerator: mp.mpf, denominator: mp.mpf) -> str:
    if denominator == 0:
        return "exact_elimination" if numerator != 0 else "both_zero"
    return str(decimal(numerator / denominator, 30))


def organize(records: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for record in records:
        record_type = record["type"]
        if record_type == "chain":
            chain_id = int(record["chain"])
            require(chain_id not in chains, f"duplicate chain {chain_id}")
            chains[chain_id] = {
                "chain": record,
                "levels": {},
                "q_terms": {},
                "steps": {},
                "end": None,
            }
        elif record_type in {"level", "q_terms", "recurrence", "chain_end"}:
            chain_id = int(record["chain"])
            require(chain_id in chains, f"record before chain header: {chain_id}")
            if record_type == "level":
                chains[chain_id]["levels"][int(record["level"])] = record
            elif record_type == "q_terms":
                chains[chain_id]["q_terms"][int(record["nit"])] = record
            elif record_type == "recurrence":
                chains[chain_id]["steps"][int(record["nit"])] = record
            else:
                chains[chain_id]["end"] = record
    require(chains, "no telemetry chains")
    return chains


def analyze_chain(chain_id: int, data: dict[str, Any]) -> dict[str, Any]:
    header = data["chain"]
    levels: dict[int, dict[str, Any]] = data["levels"]
    q_terms: dict[int, dict[str, Any]] = data["q_terms"]
    steps: dict[int, dict[str, Any]] = data["steps"]
    mit = int(header["mit"])
    require(set(levels) == set(range(1, mit + 1)), f"chain {chain_id} level gap")
    require(set(steps) == set(range(1, mit)), f"chain {chain_id} recurrence gap")
    require(set(q_terms) == set(steps), f"chain {chain_id} q/recurrence gap")
    require(data["end"] is not None, f"chain {chain_id} end missing")

    tpm = mp.mpf(header["tpm"])
    direct_ladder: dict[int, dict[int, mp.mpc]] = {}
    for level_id, level in levels.items():
        direct_ladder[level_id] = {
            dps: independent_direct_sum(level, tpm, dps) for dps in PRECISION_LADDER
        }

    direct = {
        level_id: adapt(levels[level_id], ladder[PRECISION_LADDER[-1]])
        for level_id, ladder in direct_ladder.items()
    }
    precision_drifts = [
        norm(ladder[PRECISION_LADDER[-1]] - ladder[PRECISION_LADDER[-2]])
        for ladder in direct_ladder.values()
    ]

    source_step_errors: list[mp.mpf] = []
    logged_model_errors: list[mp.mpf] = []
    logged_direct_errors: list[mp.mpf] = []
    local_defects: dict[int, mp.mpc] = {}
    logged_local_defect_errors: list[mp.mpf] = []
    q_identity_errors: list[mp.mpf] = []
    w1_norms: list[mp.mpf] = []
    q_component_sums: list[mp.mpf] = []
    q_norms: list[mp.mpf] = []

    for nit in range(mit - 1, 0, -1):
        step = steps[nit]
        q_record = q_terms[nit]
        source_before = mpc(step["state_before"])
        source_after = mpc(step["state_after"])
        source_step_errors.append(norm(source_after - apply_step(step, source_before)))

        model = apply_step(step, direct[nit + 1])
        defect = direct[nit] - model
        local_defects[nit] = defect
        logged_model_errors.append(norm(model - mpc(step["model_after"])))
        logged_local_defect_errors.append(norm(defect - mpc(step["local_defect"])))
        logged_direct_errors.extend(
            [
                norm(direct[nit] - mpc(step["adapted_parent"])),
                norm(direct[nit + 1] - mpc(step["adapted_child"])),
            ]
        )

        t1 = mpc(q_record["t1"])
        t2 = mpc(q_record["t2"])
        t3 = mpc(q_record["t3"])
        t4 = mpc(q_record["t4"])
        t5 = mpc(q_record["t5"])
        qq = mpc(q_record["qq"])
        reconstructed_q = mp.conj(t1 + t2 + t4) + t3 + t5
        q_identity_errors.append(norm(qq - reconstructed_q))
        w1_norms.append(norm(t5))
        q_component_sums.append(sum(norm(value) for value in (t1, t2, t3, t4, t5)))
        q_norms.append(norm(qq))

    final_source = mpc(data["end"]["final_state"])
    root_error = direct[1] - final_source

    if mit == 1:
        source_kernel = final_source
    else:
        source_kernel = mpc(steps[mit - 1]["state_before"])
    propagated_error = direct[mit] - source_kernel
    kernel_adapter_error = propagated_error
    for nit in range(mit - 1, 0, -1):
        propagated_error = apply_homogeneous(steps[nit], propagated_error) + local_defects[nit]
    accumulation_identity_error = norm(root_error - propagated_error)

    shell_frontier = []
    for shell_level in range(mit, 0, -1):
        shell_state = direct[shell_level]
        for nit in range(shell_level - 1, 0, -1):
            shell_state = apply_step(steps[nit], shell_state)
        residual = direct[1] - shell_state
        shell_frontier.append(
            {
                "shell_level": shell_level,
                "direct_terms": int(levels[shell_level]["length"]) + 1,
                "root_residual_norm": str(decimal(norm(residual), 30)),
                "source_error_reduction_factor": reduction_text(
                    norm(root_error), norm(residual)
                ),
                "bypassed_local_defects": list(range(shell_level, mit)),
                "bypasses_kernel_adapter": True,
            }
        )

    first_parent = next(
        (row for row in shell_frontier if row["shell_level"] == mit - 1), None
    )
    first_parent_reduction = (
        None
        if first_parent is None
        else reduction_text(
            norm(root_error), mp.mpf(first_parent["root_residual_norm"])
        )
    )
    cancellation_ratios = [
        finite_ratio(component_sum, q_norm_value)
        for component_sum, q_norm_value in zip(q_component_sums, q_norms)
    ]
    cancellation_ratios = [value for value in cancellation_ratios if value is not None]

    return {
        "chain": chain_id,
        "block": int(header["block"]),
        "sum_index": int(header["sum_index"]),
        "branch": int(header["branch"]),
        "mit": mit,
        "ip": int(header["ip"]),
        "initial_length": int(header["initial_length"]),
        "kernel_length": int(header["kernel_length"]),
        "level_lengths": [int(levels[level]["length"]) for level in range(1, mit + 1)],
        "level_degrees": [int(levels[level]["degree"]) for level in range(1, mit + 1)],
        "max_precision_ladder_drift": str(decimal(max(precision_drifts, default=mp.mpf(0)), 30)),
        "max_source_step_reconstruction_error": str(decimal(max(source_step_errors, default=mp.mpf(0)), 30)),
        "max_logged_model_reconstruction_error": str(decimal(max(logged_model_errors, default=mp.mpf(0)), 30)),
        "max_logged_direct_reconstruction_error": str(decimal(max(logged_direct_errors, default=mp.mpf(0)), 30)),
        "max_logged_local_defect_reconstruction_error": str(decimal(max(logged_local_defect_errors, default=mp.mpf(0)), 30)),
        "max_q_identity_error": str(decimal(max(q_identity_errors, default=mp.mpf(0)), 30)),
        "kernel_adapter_error_norm": str(decimal(norm(kernel_adapter_error), 30)),
        "source_root_error_norm": str(decimal(norm(root_error), 30)),
        "accumulation_identity_error": str(decimal(accumulation_identity_error, 30)),
        "local_defect_norms": {
            str(nit): str(decimal(norm(value), 30))
            for nit, value in sorted(local_defects.items())
        },
        "max_local_defect_norm": str(
            decimal(max((norm(value) for value in local_defects.values()), default=mp.mpf(0)), 30)
        ),
        "max_w1_t5_norm": str(decimal(max(w1_norms, default=mp.mpf(0)), 30)),
        "max_q_component_cancellation_ratio": str(
            decimal(max(cancellation_ratios, default=mp.mpf(0)), 30)
        ),
        "first_parent_shell": first_parent,
        "first_parent_reduction_factor": first_parent_reduction,
        "shell_frontier": shell_frontier,
    }


def aggregate(chains: list[dict[str, Any]]) -> dict[str, Any]:
    local_defects = [
        mp.mpf(value)
        for chain in chains
        for value in chain["local_defect_norms"].values()
    ]
    root_errors = [mp.mpf(chain["source_root_error_norm"]) for chain in chains]
    kernel_errors = [mp.mpf(chain["kernel_adapter_error_norm"]) for chain in chains]
    w1_norms = [mp.mpf(chain["max_w1_t5_norm"]) for chain in chains]
    first_parent_costs = [
        chain["first_parent_shell"]["direct_terms"]
        for chain in chains
        if chain["first_parent_shell"] is not None
    ]
    first_parent_residuals = [
        mp.mpf(chain["first_parent_shell"]["root_residual_norm"])
        for chain in chains
        if chain["first_parent_shell"] is not None
    ]
    first_parent_reductions = [
        mp.mpf(chain["first_parent_reduction_factor"])
        for chain in chains
        if chain["first_parent_reduction_factor"] not in {None, "exact_elimination", "both_zero"}
    ]
    first_parent_exact_eliminations = sum(
        chain["first_parent_reduction_factor"] == "exact_elimination" for chain in chains
    )

    def summary(values: list[mp.mpf]) -> dict[str, str]:
        if not values:
            return {"min": "0", "median": "0", "p90": "0", "max": "0"}
        return {
            "min": str(decimal(min(values), 30)),
            "median": str(decimal(quantile(values, 0.5), 30)),
            "p90": str(decimal(quantile(values, 0.9), 30)),
            "max": str(decimal(max(values), 30)),
        }

    def summarize_chains(selected: list[dict[str, Any]]) -> dict[str, Any]:
        defects = [mp.mpf(chain["max_local_defect_norm"]) for chain in selected]
        roots = [mp.mpf(chain["source_root_error_norm"]) for chain in selected]
        kernels = [mp.mpf(chain["kernel_adapter_error_norm"]) for chain in selected]
        w1_values = [mp.mpf(chain["max_w1_t5_norm"]) for chain in selected]
        cancellation = [
            mp.mpf(chain["max_q_component_cancellation_ratio"]) for chain in selected
        ]
        return {
            "chain_count": len(selected),
            "local_defect_norms": summary(defects),
            "source_root_error_norms": summary(roots),
            "kernel_adapter_error_norms": summary(kernels),
            "w1_t5_norms": summary(w1_values),
            "w1_to_local_defect_ratios": summary(
                [w1 / defect for w1, defect in zip(w1_values, defects) if defect != 0]
            ),
            "local_defect_to_root_error_ratios": summary(
                [defect / root for defect, root in zip(defects, roots) if root != 0]
            ),
            "kernel_to_root_error_ratios": summary(
                [kernel / root for kernel, root in zip(kernels, roots) if root != 0]
            ),
            "q_component_cancellation_ratios": summary(cancellation),
        }

    all_chain_summary = summarize_chains(chains)
    branches = sorted({chain["branch"] for chain in chains})

    return {
        "chain_count": len(chains),
        "recurrence_count": sum(len(chain["local_defect_norms"]) for chain in chains),
        "branches": branches,
        "initial_length_range": [
            min(chain["initial_length"] for chain in chains),
            max(chain["initial_length"] for chain in chains),
        ],
        "kernel_length_range": [
            min(chain["kernel_length"] for chain in chains),
            max(chain["kernel_length"] for chain in chains),
        ],
        "local_defect_norms": summary(local_defects),
        "source_root_error_norms": summary(root_errors),
        "kernel_adapter_error_norms": summary(kernel_errors),
        "w1_t5_max_per_chain_norms": summary(w1_norms),
        "w1_to_local_defect_ratios": all_chain_summary["w1_to_local_defect_ratios"],
        "local_defect_to_root_error_ratios": all_chain_summary["local_defect_to_root_error_ratios"],
        "kernel_to_root_error_ratios": all_chain_summary["kernel_to_root_error_ratios"],
        "q_component_cancellation_ratios": all_chain_summary["q_component_cancellation_ratios"],
        "by_branch": {
            str(branch): summarize_chains(
                [chain for chain in chains if chain["branch"] == branch]
            )
            for branch in branches
        },
        "first_parent_shell_term_range": (
            [min(first_parent_costs), max(first_parent_costs)] if first_parent_costs else []
        ),
        "first_parent_root_residual_norms": summary(first_parent_residuals),
        "first_parent_error_reduction_factors": summary(first_parent_reductions),
        "first_parent_exact_elimination_count": first_parent_exact_eliminations,
        "max_precision_ladder_drift": str(
            decimal(max(mp.mpf(chain["max_precision_ladder_drift"]) for chain in chains), 30)
        ),
        "max_accumulation_identity_error": str(
            decimal(max(mp.mpf(chain["accumulation_identity_error"]) for chain in chains), 30)
        ),
    }


def route_comparison(summary: dict[str, Any]) -> dict[str, Any]:
    return {
        "termwise_w1": {
            "observed_quantity": "Fortran t5, mapped to paper W1, on the saved finite chains",
            "diagnostic": summary["w1_t5_max_per_chain_norms"],
            "w1_to_local_defect_ratio": summary["w1_to_local_defect_ratios"],
            "strength": "Exposes which endpoint/saddle component is large and where cancellation is required.",
            "blocker": "Observed W1-to-defect ratios show substantial cancellation, while sample magnitudes provide no uniform contour, tail, special-function, or selector-margin bound.",
            "status": "diagnostic_only",
        },
        "direct_local_defect": {
            "observed_quantity": "Exact finite parent minus the recorded affine child recurrence",
            "diagnostic": summary["local_defect_norms"],
            "strength": "Bundles all q-component cancellation into one finite quantity and admits direct interval enclosure once the transformed-level adapter is formalized.",
            "blocker": "The present sums use rounded logged coefficients and point values, not interval coefficients over the physical shifted interval or disk.",
            "status": "best_next_certificate_candidate",
        },
        "exact_shell": {
            "observed_quantity": "Direct start at the first parent above the kernel, followed by unchanged upper recurrence",
            "term_range": summary["first_parent_shell_term_range"],
            "residual_diagnostic": summary["first_parent_root_residual_norms"],
            "reduction_diagnostic": summary["first_parent_error_reduction_factors"],
            "exact_elimination_count": summary["first_parent_exact_elimination_count"],
            "strength": "Exactly removes the kernel adapter and first post-kernel local defect without separately bounding them.",
            "blocker": "All saved chains have MIT=2, so the first parent is the original 105-term sum; exact elimination here is not yet evidence of a cheaper physical-height shell.",
            "status": "conditional_algorithmic_bypass_candidate",
        },
        "provisional_order": [
            "direct_local_defect",
            "exact_shell",
            "termwise_w1",
        ],
        "selection_rule": (
            "Promote the direct-defect route first if interval re-summation remains affordable; "
            "use the first-parent shell when it removes the dominant deepest error at acceptable "
            "measured cost; return to termwise W1 only where a uniform analytic majorant is needed "
            "to avoid direct summation."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate_result = artifact["aggregate"]
    routes = artifact["route_comparison"]
    return "\n".join(
        [
            "# Hardy Chain Telemetry Scout",
            "",
            "Date: 2026-08-06",
            "Status: finite low-height diagnostic; not a proof and no external error theorem",
            "",
            "## Evaluator Boundary",
            "",
            "The accepted resumable evaluator remains byte-pinned. The separately compiled telemetry derivative matches its full-precision checkpoint journal after removing only the run provenance id, and enabling telemetry leaves that derivative's journal and displayed values exactly unchanged.",
            "",
            "The phase constant is not inserted as an unexplained `pi`. The accepted source computes `p = 4*atan(1)`, then `tpm = -2*p`; each chain records that actual `tpm`, and every independent reconstruction consumes the recorded value.",
            "",
            "## Independent Reconstruction",
            "",
            f"The scout contains `{aggregate_result['chain_count']}` chains and `{aggregate_result['recurrence_count']}` recurrence rows. Every finite transformed level is re-summed independently at {', '.join(str(value) for value in PRECISION_LADDER)} decimal digits from its logged length and coefficients.",
            "",
            f"Initial lengths range from `{aggregate_result['initial_length_range'][0]}` to `{aggregate_result['initial_length_range'][1]}` terms before the inclusive endpoint adjustment; kernel lengths range from `{aggregate_result['kernel_length_range'][0]}` to `{aggregate_result['kernel_length_range'][1]}`. The maximum 80-to-120-digit ladder drift is `{aggregate_result['max_precision_ladder_drift']}`.",
            "",
            f"The conjugation-aware affine accumulation identity reconstructs the observed root error with maximum residual `{aggregate_result['max_accumulation_identity_error']}`.",
            "",
            "## Route Comparison",
            "",
            f"- **Direct local defect:** observed norm summary `{json.dumps(aggregate_result['local_defect_norms'], sort_keys=True)}`. The W1-to-defect ratio summary is `{json.dumps(aggregate_result['w1_to_local_defect_ratios'], sort_keys=True)}`. This is the leading certificate candidate because it preserves cancellations instead of bounding five source components separately.",
            f"- **Exact first-parent shell:** direct term range `{routes['exact_shell']['term_range']}` and exact elimination on `{routes['exact_shell']['exact_elimination_count']}` chains. Because every saved chain has `MIT=2`, this shell is the original 105-term sum; it validates the bypass identity but does not yet demonstrate a cheaper scalable shell.",
            f"- **Termwise W1:** observed `t5` summary `{json.dumps(aggregate_result['w1_t5_max_per_chain_norms'], sort_keys=True)}`. It remains necessary as an analytic fallback, but these samples supply no uniform W1 bound.",
            "",
            "Provisional order: direct finite-defect enclosure, first-parent exact shell, then termwise W1 where direct cost forces it.",
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    for path in (FIXTURE_RESULT, TELEMETRY):
        require(path.is_file(), f"missing telemetry fixture artifact: {path}")
    fixture = json.loads(FIXTURE_RESULT.read_text(encoding="utf-8"))
    require(fixture["enabled_vs_disabled_journal_exact"], "telemetry state equivalence missing")
    require(
        fixture["accepted_reference_journal_exact_except_run_id"],
        "accepted-reference state equivalence missing",
    )

    records = load_jsonl(TELEMETRY)
    chains = organize(records)
    analyses = [analyze_chain(chain_id, chains[chain_id]) for chain_id in sorted(chains)]
    aggregate_result = aggregate(analyses)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_chain_telemetry_scout",
        "status": "finite_low_height_route_comparison",
        "phase_constant_provenance": {
            "accepted_source_definition": "p=4*ATAN(1); tpp=2*p; tpm=-tpp",
            "reconstruction_input": "Each chain's recorded tpm; no independently inserted pi constant.",
        },
        "precision_ladder_decimal_digits": list(PRECISION_LADDER),
        "evaluator_equivalence": {
            "enabled_vs_disabled_full_precision_journal_exact": True,
            "accepted_reference_journal_exact_except_run_id": True,
            "enabled_vs_disabled_displayed_output_exact": True,
            "accepted_reference_displayed_output_exact": True,
            "fixture_path": relative(FIXTURE_RESULT),
            "fixture_sha256": file_hash(FIXTURE_RESULT),
        },
        "aggregate": aggregate_result,
        "route_comparison": route_comparison(aggregate_result),
        "chains": analyses,
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "builder": {"path": relative(SCRIPT), "sha256": file_hash(SCRIPT)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_executable_target": (
            "Replace point coefficients by outward-rounded coefficient intervals for the "
            "smallest-cost direct-defect and first-parent-shell rows; rerun over the entire "
            "15-shift fixture interval, then inventory selector and hierarchy transition margins."
        ),
        "proof_boundary": (
            "This scout proves evaluator non-interference for one accepted low-height fixture, "
            "independently reconstructs its finite transformed sums at an increasing precision "
            "ladder, verifies the conjugation-aware finite defect accumulation numerically, and "
            "measures exact-shell costs. It does not prove a transformed-level interval adapter, "
            "uniform W1 or special-function bounds, selector/transition exclusion, outer Hardy "
            "representation error, physical-height complexity, evaluator C6 or disk control, a "
            "carrier value, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }

    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    result_text = json.dumps(artifact, indent=2) + "\n"
    note_text = render_note(artifact)
    result_tmp = RESULT.with_suffix(".json.tmp")
    note_tmp = NOTE.with_suffix(".md.tmp")
    result_tmp.write_text(result_text, encoding="utf-8")
    note_tmp.write_text(note_text, encoding="utf-8")
    result_tmp.replace(RESULT)
    note_tmp.replace(NOTE)
    print(
        "analyzed Hardy chain telemetry routes: "
        f"{aggregate_result['chain_count']} chains, "
        f"{aggregate_result['recurrence_count']} recurrence rows, "
        "three diagnostic routes compared"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
