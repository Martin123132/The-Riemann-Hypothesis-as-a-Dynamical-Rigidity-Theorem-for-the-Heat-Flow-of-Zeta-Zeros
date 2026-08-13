#!/usr/bin/env python3
"""Prove the exact reflection pairing of the two 399-mode branch rosters."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_opposite_branch_roster_reflection_pairing_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "event_selection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_hankel_carrier_event_selection_gate.json",
    "branch_envelope": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate.json",
    "event_atlas": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.json",
    "selector_reassembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate.json",
    "boundary_fourier_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.json",
}

OLD_C = 159_577
NEW_C = 159_579
Q = 39_894
OLD_LO = 39_695
OLD_HI = 40_093
NEW_LO = 39_696
NEW_HI = 40_094
PAIR_COUNT = 199
PRECISION = 100


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def roster_certificate() -> dict[str, Any]:
    old_modes = list(range(OLD_LO, OLD_HI + 1))
    new_modes = list(range(NEW_LO, NEW_HI + 1))
    old_labels = [4 * mode - OLD_C for mode in old_modes]
    new_labels = [4 * mode - NEW_C for mode in new_modes]
    boundary_labels = [4 * mode - OLD_C for mode in new_modes]

    require(len(old_modes) == len(new_modes) == 399, "selector roster length drift")
    require(old_labels == list(range(-797, 796, 4)), "old detuning roster drift")
    require(new_labels == list(range(-795, 798, 4)), "adjacent detuning roster drift")
    require(boundary_labels == list(range(-793, 800, 4)), "boundary Fourier detuning roster drift")

    old_pairs = [(-(4 * j + 1), 4 * j + 3) for j in range(PAIR_COUNT)]
    new_pairs = [(-(4 * j + 3), 4 * j + 1) for j in range(PAIR_COUNT)]
    boundary_pairs = [(-(4 * j + 1), 4 * j + 3) for j in range(PAIR_COUNT)]
    require({value for pair in old_pairs for value in pair} | {-797} == set(old_labels), "old reflection pairing is not exhaustive")
    require({value for pair in new_pairs for value in pair} | {797} == set(new_labels), "adjacent reflection pairing is not exhaustive")
    require({value for pair in boundary_pairs for value in pair} | {799} == set(boundary_labels), "boundary reflection pairing is not exhaustive")

    old_only = sorted(set(old_modes) - set(new_modes))
    new_only = sorted(set(new_modes) - set(old_modes))
    require(old_only == [OLD_LO] and new_only == [NEW_HI], "selector edge exchange drift")
    require(4 * OLD_LO - OLD_C == -797, "old unpaired edge label drift")
    require(4 * NEW_HI - OLD_C == 799, "old-coordinate missing mate drift")
    require(4 * NEW_HI - NEW_C == 797, "adjacent unpaired edge label drift")

    return {
        "old_selector": {
            "C": OLD_C,
            "mode_roster": [OLD_LO, OLD_HI],
            "detuning_label_roster": [-797, 795, 4],
            "negative_count": 200,
            "positive_count": 199,
            "reflection_pair_count": PAIR_COUNT,
            "unpaired_mode": OLD_LO,
            "unpaired_label": -797,
            "missing_positive_mate_mode": NEW_HI,
            "missing_positive_mate_label_in_old_coordinates": 799,
        },
        "adjacent_selector": {
            "C": NEW_C,
            "mode_roster": [NEW_LO, NEW_HI],
            "detuning_label_roster": [-795, 797, 4],
            "negative_count": 199,
            "positive_count": 200,
            "reflection_pair_count": PAIR_COUNT,
            "unpaired_mode": NEW_HI,
            "unpaired_label": 797,
            "missing_negative_mate_mode": OLD_LO,
            "missing_negative_mate_label_in_adjacent_coordinates": -799,
        },
        "shared_boundary_fourier_block": {
            "strip_lower_endpoint_C": OLD_C,
            "mode_roster": [NEW_LO, NEW_HI],
            "detuning_label_roster_in_old_strip_coordinates": [-793, 799, 4],
            "negative_count": 199,
            "positive_count": 200,
            "reflection_pair_count": PAIR_COUNT,
            "unpaired_mode": NEW_HI,
            "unpaired_label": 799,
            "coordinate_warning": "Same mode set as the adjacent lower-edge event roster, but detunings use the old one-cell strip scale and labels 4m-C, not 4m-(C+2).",
        },
        "old_lower_to_boundary_mode_set_exchange": {
            "common_mode_roster": [NEW_LO, OLD_HI],
            "dropped_mode": OLD_LO,
            "added_mode": NEW_HI,
            "scope": "Finite mode-set identity only; not the complete selector jump or a chart-amplitude identity.",
        },
    }


def symbolic_certificate() -> dict[str, str]:
    q = sp.symbols("q", nonzero=True)
    c_plus, c_minus = sp.symbols("C_plus C_minus")

    old_direct = c_plus * sum(q ** (4 * j + 1) for j in range(200))
    old_direct += c_minus * sum(q ** (-(4 * j + 3)) for j in range(199))
    old_paired = c_plus * q**797
    old_paired += sum(
        c_plus * q ** (4 * j + 1) + c_minus * q ** (-(4 * j + 3))
        for j in range(PAIR_COUNT)
    )
    require(sp.expand(old_direct - old_paired) == 0, "old branch pairing identity failed")

    new_direct = c_plus * sum(q ** (4 * j + 3) for j in range(199))
    new_direct += c_minus * sum(q ** (-(4 * j + 1)) for j in range(200))
    new_paired = c_minus * q**-797
    new_paired += sum(
        c_plus * q ** (4 * j + 3) + c_minus * q ** (-(4 * j + 1))
        for j in range(PAIR_COUNT)
    )
    require(sp.expand(new_direct - new_paired) == 0, "adjacent branch pairing identity failed")

    old_negative_kernel = q * (1 - q**800) / (1 - q**4)
    old_positive_kernel = q**-3 * (1 - q**-796) / (1 - q**-4)
    require(sp.cancel(old_direct - (c_plus * old_negative_kernel + c_minus * old_positive_kernel)) == 0, "old half-kernel factorization failed")

    new_negative_kernel = q**3 * (1 - q**796) / (1 - q**4)
    new_positive_kernel = q**-1 * (1 - q**-800) / (1 - q**-4)
    require(sp.cancel(new_direct - (c_plus * new_negative_kernel + c_minus * new_positive_kernel)) == 0, "adjacent half-kernel factorization failed")

    boundary_direct = c_plus * sum(q ** (4 * j + 1) for j in range(199))
    boundary_direct += c_minus * sum(q ** (-(4 * j + 3)) for j in range(200))
    boundary_paired = c_minus * q**-799
    boundary_paired += sum(
        c_plus * q ** (4 * j + 1) + c_minus * q ** (-(4 * j + 3))
        for j in range(PAIR_COUNT)
    )
    require(sp.expand(boundary_direct - boundary_paired) == 0, "boundary branch pairing identity failed")
    boundary_negative_kernel = q * (1 - q**796) / (1 - q**4)
    boundary_positive_kernel = q**-3 * (1 - q**-800) / (1 - q**-4)
    require(sp.cancel(boundary_direct - (c_plus * boundary_negative_kernel + c_minus * boundary_positive_kernel)) == 0, "boundary half-kernel factorization failed")

    return {
        "branch_definition": "C_+(y)=exp(i*(xi-pi/4))*S_+(X,y); C_-(y)=conj(C_+(y)) for real X,U,V",
        "old_opposite_branch_sum": "O_C=C_+*sum_(j=0)^199 exp(i*(4j+1)*a*y)+C_-*sum_(j=0)^198 exp(-i*(4j+3)*a*y)",
        "old_half_kernel_form": "O_C=C_+*exp(i*399*x)*sin(400*x)/sin(2*x)+C_-*exp(-i*399*x)*sin(398*x)/sin(2*x), x=a*y",
        "old_reflection_form": "O_C=exp(i*797*x)C_+ + exp(-i*x)*sin(398*x)/sin(2*x)*(exp(i*398*x)C_+ + exp(-i*398*x)C_-)",
        "old_real_projection_form": "O_C=exp(i*797*x)C_+ + 2*exp(-i*x)*sin(398*x)/sin(2*x)*Re(exp(i*398*x)C_+)",
        "adjacent_opposite_branch_sum": "O_(C+2)=C_+'*sum_(j=0)^198 exp(i*(4j+3)*a'*y)+C_-'*sum_(j=0)^199 exp(-i*(4j+1)*a'*y)",
        "adjacent_half_kernel_form": "O_(C+2)=C_+'*exp(i*399*x')*sin(398*x')/sin(2*x')+C_-'*exp(-i*399*x')*sin(400*x')/sin(2*x')",
        "adjacent_reflection_form": "O_(C+2)=exp(-i*797*x')C_-' + exp(i*x')*sin(398*x')/sin(2*x')*(exp(i*398*x')C_+' + exp(-i*398*x')C_-')",
        "removable_zero_rule": "Every sine quotient is defined by its finite exponential sum at sin(2*x)=0.",
        "old_endpoint_value": "O_C(0)=200*C_+ + 199*C_-=399*Re(C_+)+i*Im(C_+)",
        "adjacent_endpoint_value": "O_(C+2)(0)=199*C_+'+200*C_-'=399*Re(C_+')-i*Im(C_+')",
        "boundary_opposite_branch_sum": "O_F=C_+*sum_(j=0)^198 exp(i*(4j+1)*x)+C_-*sum_(j=0)^199 exp(-i*(4j+3)*x)",
        "boundary_half_kernel_form": "O_F=C_+*exp(i*397*x)*sin(398*x)/sin(2*x)+C_-*exp(-i*401*x)*sin(400*x)/sin(2*x)",
        "boundary_reflection_form": "O_F=exp(-i*799*x)C_- + exp(-i*x)*sin(398*x)/sin(2*x)*(exp(i*398*x)C_+ + exp(-i*398*x)C_-)",
        "boundary_real_projection_form": "O_F=exp(-i*799*x)C_- + 2*exp(-i*x)*sin(398*x)/sin(2*x)*Re(exp(i*398*x)C_+)",
        "boundary_endpoint_value": "O_F(0)=199*C_+ + 200*C_-=399*Re(C_+)-i*Im(C_+)",
    }


def scale_certificate() -> dict[str, str]:
    ctx.dps = PRECISION
    pi = arb.pi()
    beta_old = (pi * arb(OLD_C) ** 2 / 8) ** (arb(1) / 3)
    beta_new = (pi * arb(NEW_C) ** 2 / 8) ** (arb(1) / 3)
    a_old = beta_old / OLD_C
    a_new = beta_new / NEW_C
    mismatch = abs(a_new - a_old)
    require(mismatch < arb("6e-8"), "adjacent detuning scale mismatch ceiling failed")
    require(mismatch > 0, "adjacent detuning scales unexpectedly identical")
    return {
        "old_beta_ball": beta_old.str(PRECISION, more=True),
        "adjacent_beta_ball": beta_new.str(PRECISION, more=True),
        "old_a_equals_beta_over_C_ball": a_old.str(PRECISION, more=True),
        "adjacent_a_equals_beta_over_C_ball": a_new.str(PRECISION, more=True),
        "absolute_a_scale_mismatch_ball": mismatch.str(PRECISION, more=True),
        "scale_mismatch_below": "6e-8",
    }


def render_note(artifact: dict[str, Any]) -> str:
    s = artifact["symbolic_certificate"]
    scale = artifact["scale_certificate"]
    return f"""# Beta^-4 opposite-branch roster reflection pairing

Date: 2026-08-13
Status: exact finite-roster symmetry; not a proof of the grouped integral

Put `a=beta/C`, `x=a*y`, and write the exact branches from Sections
11.378--11.380 as

```text
C_+(y)=exp(i[xi-pi/4])S_+(X,y),
C_-(y)=conj(C_+(y)).                                  (RP1)
```

For the fixed selector `C={OLD_C}`, the 399 modes `{OLD_LO}..{OLD_HI}` have
detuning labels `4m-C=-797,-793,...,795`.  There are 200 negative modes and
199 positive modes.  The branch opposite to the endpoint crossing is `C_+`
on the negative half and `C_-` on the positive half, so its complete roster
sum is

```text
{s['old_opposite_branch_sum']}.                        (RP2)
```

Pair the labels `-(4j+1)` and `+(4j+3)` for `0<=j<=198`.  The only unpaired
old-selector label is `-797`, belonging to mode `{OLD_LO}`.  Exact geometric
summation gives

```text
{s['old_half_kernel_form']},                           (RP3)

{s['old_real_projection_form']}.                      (RP4)
```

All quotients in (RP3)--(RP4) are removable: at `sin(2x)=0` their values are
defined by the finite sums in (RP2).  Thus (RP4) is 199 conjugate pairs plus
one edge mode, equivalently one real Dirichlet projection plus that edge,
not 399 unrelated branch errors.

At the shared selector boundary, the already-certified one-cell Fourier block
has modes `{NEW_LO}..{NEW_HI}` but still uses the old strip coordinate `C`.
Its detuning labels are therefore `-793,-789,...,799`, not the adjacent-chart
labels.  Its exact reflected form is

```text
{s['boundary_half_kernel_form']},                     (RP5)

{s['boundary_reflection_form']}.                     (RP6)
```

This block has 199 negative modes, 200 positive modes, and unpaired label
`+799`, mode `{NEW_HI}`.  Set-theoretically, passing from the old lower-edge
roster to this boundary block keeps `{NEW_LO}..{OLD_HI}`, drops `{OLD_LO}`,
and adds `{NEW_HI}`.  This is only a finite mode-set identity; the complete
selector jump also contains the symmetric Poisson complement and endpoint
half-current already certified by the Fourier-completion gate.

The adjacent selector `C+2={NEW_C}` has the same mode set but its own event
coordinate, with labels `-795,-791,...,797`.  Its separate mirror formula is

```text
{s['adjacent_half_kernel_form']},                      (RP7)

{s['adjacent_reflection_form']}.                      (RP8)
```

The adjacent unpaired label is `+797`, again mode `{NEW_HI}`.  At `y=0`, the
old lower-edge and shared-boundary projections have opposite one-copy
imaginary imbalances,

```text
{s['old_endpoint_value']},
{s['boundary_endpoint_value']}.                       (RP9)
```

The two canonical charts do not have identical scales.  Direct interval
evaluation gives

```text
a_C={scale['old_a_equals_beta_over_C_ball']},
a_(C+2)={scale['adjacent_a_equals_beta_over_C_ball']},
|a_(C+2)-a_C|={scale['absolute_a_scale_mismatch_ball']}<6e-8. (RP10)
```

The boundary Fourier-completion theorem already combines its whole 399 block,
the zero/negative/outer-positive complement, and the endpoint half-current
before absolute values.  Therefore (RP6) is an internal branch organization,
not a license to replace that completed identity by an edge difference.  The
next admissible use is to insert (RP6) into the completed object and determine
which branch projection continues into the ordinary-Morse corridor.  A
termwise absolute sum over 399 modes would erase this exact structure.

This gate proves the three roster reflection algebras only.  It does not bound
the grouped opposite-branch integral, identify the old and adjacent chart
amplitudes, replace the completed Poisson identity by edge cancellation,
close `Q_K-T` or `T_upper`.  No height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    dependencies: dict[str, Any] = {}
    for name, path in DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {name}")
        dependencies[name] = json.loads(path.read_text(encoding="utf-8"))

    require(dependencies["event_selection"]["decision"]["opposite_extracted_carrier_nonstationary_on_fold_strip"] is True, "event selection drift")
    require(dependencies["branch_envelope"]["decision"]["exact_scaled_branch_amplitude_bound_proved"] is True, "branch envelope drift")
    atlas = dependencies["event_atlas"]["certified_event_atlas"]
    require(atlas["lower_edge_transition_roster"] == [OLD_LO, OLD_HI], "old atlas roster drift")
    require(atlas["adjacent_selector_lower_edge_transition_roster"] == [NEW_LO, NEW_HI], "adjacent atlas roster drift")
    require(dependencies["selector_reassembly"]["decision"]["complete_symmetric_poisson_strip_reassembles_removed_term"] is True, "selector reassembly drift")
    require(dependencies["boundary_fourier_completion"]["decision"]["canonical_399_mode_integral_rigorously_computed"] is True, "boundary Fourier completion drift")
    require(dependencies["boundary_fourier_completion"]["geometry"]["stationary_mode_range"] == [NEW_LO, NEW_HI], "boundary Fourier roster drift")

    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_old_boundary_and_adjacent_399_mode_opposite_branch_reflection_pairings_proved",
        "passed": True,
        "roster_certificate": roster_certificate(),
        "symbolic_certificate": symbolic_certificate(),
        "scale_certificate": scale_certificate(),
        "decision": {
            "old_399_mode_opposite_branch_reduced_to_199_reflection_pairs_plus_one_edge": True,
            "adjacent_399_mode_opposite_branch_is_exact_mirror_pairing": True,
            "shared_boundary_399_mode_opposite_branch_is_exact_mirror_pairing": True,
            "old_lower_to_boundary_mode_set_exchanges_two_edges": True,
            "edge_exchange_alone_reassembles_selector_jump": False,
            "termwise_399_mode_absolute_sum_required": False,
            "grouped_opposite_branch_integral_bound_proved": False,
            "old_to_adjacent_chart_amplitude_identification_proved": False,
            "edge_plus_endpoint_half_current_cancellation_proved": False,
            "complete_Q_K_minus_T_bound_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "next_action": "Insert the shared-boundary reflection formula into the already-completed Fourier object and derive which grouped branch projection continues into the ordinary-Morse corridor without separating the Poisson complement or endpoint half-current.",
        "proof_boundary": "Exact finite-roster reflection pairings and a finite mode-set exchange only. No grouped integral bound, chart-amplitude match, replacement of the completed Poisson identity by edge cancellation, Q_K-T bound, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
        "dependencies": {name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "resource_policy": {"workers": 1, "process_priority": priority},
        "runtime": {"elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    artifact["sources"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("proved 399-mode branch reflection pairing: 199 pairs plus one exchanged selector edge", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
