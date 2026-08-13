#!/usr/bin/env python3
"""Independently check the source transition/endpoint scheduling gate."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, ctx
import mpmath as mp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.endpoint_gate.selector_gate.event_gate.window_gate.cell_gate.ode_gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing transition scheduling artifact")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "transition scheduling artifact does not pass")

    for key, path in (
        ("upper_endpoint_gate", gate.ENDPOINT_GATE),
        ("selector_strip_gate", gate.SELECTOR_GATE),
        ("pinned_source", gate.SOURCE),
        ("paper_2015", gate.PAPER_2015),
        ("paper_2026", gate.PAPER_2026),
    ):
        require(artifact["dependencies"][key]["sha256"] == gate.file_hash(path), f"{key} hash mismatch")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")

    require(artifact["source_contract"] == gate.source_contract(), "source contract mismatch")
    ctx.dps = 130
    fresh = gate.event_interval_certificate()
    require(artifact["certificate"] == fresh, "interval certificate mismatch")

    c = artifact["certificate"]
    pre = list(range(c["pre_effective_start"], c["pre_endpoint"] + 1, 2))
    post = list(range(c["post_effective_start"], c["post_endpoint"] + 1, 2))
    fixed = list(range(c["fixed_chart_start"], c["fixed_chart_endpoint"] + 1, 2))
    require(len(pre) == len(post) == c["pre_roster_count"], "independent roster count failed")
    require(len(fixed) == len(pre) + 1, "fixed-chart count failed")
    require(set(pre) - set(post) == {gate.C}, "lower roster difference failed")
    require(set(post) - set(pre) == {gate.B_PLUS}, "upper roster difference failed")
    require(set(fixed) - set(pre) == {gate.B_PLUS}, "pre/fixed chart relation failed")
    require(set(fixed) - set(post) == {gate.C}, "post/fixed chart relation failed")

    # A fresh numerical implementation is only a diagnostic consistency check.
    mp.mp.dps = 100
    t64 = gate.transition_quadrature(64)
    t96 = gate.transition_quadrature(96)
    lower = gate.direct_term_point(gate.C)
    upper = gate.direct_term_point(gate.B_PLUS)
    require(abs(t64 - t96) < mp.mpf("1e-40"), "independent transition convergence failed")
    require(abs((t64 - lower) + upper) > mp.mpf("0.03"), "diagnostic source change failed")

    decisions = artifact["decision"]
    require(decisions["earlier_coupled_source_event_decomposed"] is True, "event decision corrupted")
    require(decisions["transition_is_exact_handoff_identity"] is False, "transition boundary corrupted")
    require(decisions["moving_source_endpoint_admissible_as_theorem_selector"] is False, "selector boundary corrupted")
    require(decisions["fixed_analytic_endpoint_for_next_chart"] == gate.B_PLUS, "fixed endpoint corrupted")
    require(decisions["rs_cutoff_changes_at_event"] is False, "RS cutoff decision corrupted")
    require(decisions["complete_source_error_bounded"] is False, "source-error boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated source transition scheduling independently: fixed B=5122423; priority={priority}")


if __name__ == "__main__":
    main()
