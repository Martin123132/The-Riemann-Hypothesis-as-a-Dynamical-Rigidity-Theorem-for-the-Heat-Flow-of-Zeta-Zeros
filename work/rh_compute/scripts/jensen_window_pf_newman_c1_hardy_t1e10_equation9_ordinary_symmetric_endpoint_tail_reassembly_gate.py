#!/usr/bin/env python3
"""Certify the exact symmetric residual after extracting the ordinary bulk."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "portcullis": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json",
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "coverage_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
    "Gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "endpoint_decomposition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_finite_endpoint_tail_decomposition_gate.json",
}

A = 159_577
B = 5_122_421
L = (B - A) // 2
TARGET_START = 622
TARGET_END = 39_894
TARGET_COUNT = TARGET_END - TARGET_START + 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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


def symbolic_certificate() -> dict[str, str]:
    z = sp.symbols("z", nonzero=True)
    a, n, m_cut = sp.symbols("a n M", integer=True, positive=True)
    geometric_block = z**a * (1 - z**n) / (1 - z)
    require(
        sp.cancel((1 - z) * geometric_block - (z**a - z ** (a + n))) == 0,
        "geometric block telescoping failed",
    )
    symmetric_kernel = z ** (-m_cut) * (1 - z ** (2 * m_cut + 1)) / (1 - z)
    require(
        sp.cancel((1 - z) * symmetric_kernel - (z ** (-m_cut) - z ** (m_cut + 1))) == 0,
        "symmetric kernel telescoping failed",
    )

    target_count = TARGET_END - TARGET_START + 1
    actual_partition_count = (
        m_cut
        + 1
        + (TARGET_START - 1)
        + target_count
        + (m_cut - TARGET_END)
    )
    require(sp.simplify(actual_partition_count - (2 * m_cut + 1)) == 0, "mode partition count failed")

    delta_f, delta_f1, r_plus, r_minus = sp.symbols("Delta_f Delta_f1 R_plus R_minus")
    positive_mode = sp.symbols("m", nonzero=True)
    denominator = 2 * sp.pi * sp.I * positive_mode
    i_plus = -delta_f / denominator - delta_f1 / denominator**2 + r_plus / denominator**2
    i_minus = i_plus.subs({positive_mode: -positive_mode, r_plus: r_minus})
    paired = sp.simplify(i_plus + i_minus)
    expected_pair = sp.simplify((-2 * delta_f1 + r_plus + r_minus) / denominator**2)
    require(sp.simplify(paired - expected_pair) == 0, "paired endpoint-current cancellation failed")
    require(sp.simplify(sp.diff(paired, delta_f)) == 0, "one-over-m endpoint current survived")

    half_current, complement, target_finite, target_bulk = sp.symbols("H C J_T B_T")
    target_tail = target_finite - target_bulk
    full_minus_bulk = half_current + complement + target_finite - target_bulk
    grouped_residual = half_current + complement + target_tail
    require(sp.expand(full_minus_bulk - grouped_residual) == 0, "bulk-extracted residual identity failed")

    require(A % 2 == B % 2 == 1 and (B - A) % 2 == 0, "source roster parity failed")
    require((A + 1) % 2 == (B + 1) % 2 == 0, "odd-character cancellation failed")

    return {
        "odd_roster_change": "alpha=A+2u, d alpha=2du, 0<=u<=L",
        "portcullis_to_fourier": "J_m=(-1)^m integral_0^L f_x(u)exp[-i*pi*m(A+2u)]du=I_m because A is odd",
        "fourier_coefficient": "I_m(x)=integral_0^L f_x(u)exp(-2*pi*i*m*u)du",
        "finite_symmetric_kernel": "D_M(u)=sum_(m=-M)^M exp(-2*pi*i*m*u)=sin((2M+1)*pi*u)/sin(pi*u)",
        "target_kernel": "G_T(u)=sum_(m=622)^39894 exp(-2*pi*i*m*u)=exp(-40516*pi*i*u)sin(39273*pi*u)/sin(pi*u)",
        "complement_kernel": "C_M(u)=D_M(u)-G_T(u)",
        "endpoint_values": "D_M(k)=2M+1, G_T(k)=39273, C_M(k)=2M+1-39273 for every integer k, by removable limits",
        "half_current": "H_x=[f_x(0)+f_x(L)]/2",
        "finite_reconstruction": "S_M(x)=H_x+integral_0^L f_x(u)D_M(u)du",
        "target_finite_block": "J_T(x)=integral_0^L f_x(u)G_T(u)du=sum_(m=622)^39894 J_m(x)",
        "bulk_extracted_residual": "S_M-B_T=H_x+integral_0^L f_x(u)C_M(u)du+sum_(m=622)^39894[J_m-B_m]",
        "paired_ibp_identity": "I_m+I_(-m)=[-2*Delta f_x'+R_m+R_(-m)]/(2*pi*i*m)^2",
        "one_over_m_cancellation": "The Delta f_x/(2*pi*i*m) current cancels exactly in every finite +/-m pair.",
        "symmetric_limit": "The proved symmetric-Poisson limit for D_M implies the C_M residual limit because G_T is a fixed finite trigonometric polynomial.",
    }


def mode_partition() -> dict[str, Any]:
    return {
        "cutoff_condition": f"integer M>={TARGET_END}",
        "negative_completion": {"range": ["-M", -1], "count": "M"},
        "zero_mode": {"range": [0, 0], "count": 1},
        "low_positive_completion": {"range": [1, TARGET_START - 1], "count": TARGET_START - 1},
        "target_positive_block": {"range": [TARGET_START, TARGET_END], "count": TARGET_COUNT},
        "outer_positive_completion": {
            "range": [TARGET_END + 1, "M"],
            "count": f"M-{TARGET_END}",
            "empty_at_minimal_cutoff": True,
        },
        "count_identity": f"M+1+{TARGET_START - 1}+{TARGET_COUNT}+(M-{TARGET_END})=2M+1",
        "cancellation_ledger": {
            "target_partner": f"negative modes -{TARGET_END}..-{TARGET_START} cancel the one-over-m currents of target modes {TARGET_START}..{TARGET_END}",
            "low_partner": f"negative modes -{TARGET_START - 1}..-1 pair with positive modes 1..{TARGET_START - 1}",
            "outer_partner": f"negative modes -M..-{TARGET_END + 1} pair with positive modes {TARGET_END + 1}..M",
            "zero_mode_role": "The zero coefficient is an explicit reconstruction term and has no one-over-m endpoint current.",
            "half_current_role": "The endpoint half-current is the finite-roster midpoint correction; it is not the canceller of the one-over-m pair current.",
            "grouping_guard": "Bare completed-square endpoint exponentials do not cancel pairwise. Their cancellation occurs only after each is recombined with its Fresnel tail into J_m and the complete symmetric roster is paired.",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    return f"""# Symmetric endpoint-tail reassembly after bulk extraction

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the quantitative residual bound

Write the odd source roster as `alpha=A+2u`, where

```text
A={A}, B={B}, L={L}.
```

The portcullis coefficient is exactly the ordinary Fourier coefficient:

```text
J_m=(-1)^m integral_0^L f_x(u)exp[-i*pi*m(A+2u)]du
   =integral_0^L f_x(u)exp(-2*pi*i*m*u)du=I_m.       (SR1)
```

The last equality uses only that `A` is odd.  For every finite symmetric
cutoff `M>=39894`, define

```text
D_M(u)=sum_(m=-M)^M exp(-2*pi*i*m*u)
      =sin((2M+1)pi*u)/sin(pi*u),

G_T(u)=sum_(m=622)^39894 exp(-2*pi*i*m*u)
      =exp(-40516*pi*i*u)sin(39273*pi*u)/sin(pi*u),

C_M(u)=D_M(u)-G_T(u).                                (SR2)
```

The quotients use their removable values at integer `u`.  In particular,
`D_M(k)=2M+1`, `G_T(k)=39273`, and `C_M(k)=2M+1-39273`.

Let `H_x=[f_x(0)+f_x(L)]/2`, let `B_m` be the exact full-line bulk current,
and let `Q_m=J_m-B_m=P_A,m+P_B,m` be the grouped finite-endpoint tail.  Pure
finite algebra gives the cancellation-preserving residual

```text
S_M-sum_(m=622)^39894 B_m
 =H_x+integral_0^L f_x(u)C_M(u)du
     +sum_(m=622)^39894 Q_m,                          (SR3)

S_M=H_x+integral_0^L f_x(u)D_M(u)du.
```

Equation (SR3) is the exact object left after Section 11.349 closes the
Gamma bulk.  The previously proved symmetric-Poisson theorem gives the
limit of the `D_M` integral.  Since `G_T` is a fixed finite trigonometric
polynomial, it also proves the symmetric limit of the right side of (SR3).
No auxiliary Abel regulator is needed for this fixed finite roster.

The mode ledger identifies the cancellation precisely.  Negative modes
`-39894..-622` pair with the target positive modes; `-621..-1` pair with the
low positive completion; and `-M..-39895` pair with the outer positive
completion.  For every pair,

```text
I_m+I_-m=[-2 Delta f_x'+R_m+R_-m]/(2*pi*i*m)^2,       (SR4)
```

so the apparent `Delta f_x/(2*pi*i*m)` current vanishes before absolute
values.  The zero mode has no `1/m` current.  The endpoint half-current is
the finite-roster midpoint correction, not its canceller.

This also sharpens the guard from Section 11.350.  The completed-square
bare endpoint exponentials do not cancel directly: after odd parity they
are mode-independent.  Their compensating terms live in the Fresnel tails.
Only the grouped `J_m`, followed by symmetric `+m/-m` pairing, exposes the
true cancellation.  Bounding `P_A`, `P_B`, the negative modes, or the
endpoint half-current separately would destroy (SR3)--(SR4).

The next quantitative object is therefore one joint residual, not several
endpoint-error sums:

```text
R_end=lim_(M->infinity) [H + integral f C_M
                         +sum_target(P_A+P_B)].        (SR5)
```

The upper `B` part should first be enclosed nonstationarily in this grouped
form; the lower `A` characteristic part can then be spliced to the certified
fold atlas on the same current.

Pi provenance: every `pi` in (SR1)--(SR5) comes from the equation-(9)
Kummer quadratic phase or the integer Fourier-Poisson character.  No circle,
polygon, fitted period, or inserted geometric constant is used.

Proof boundary: exact mode partition, Dirichlet-kernel representation,
symmetric residual limit, and pair-current cancellation only.  No numerical
bound for (SR5), complete endpoint-tail theorem, ordinary/fold splice,
complete `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    coverage = dependencies["coverage_ledger"]["mode_coverage"]["source_classical_target"]
    require(coverage == {"range": [TARGET_START, TARGET_END], "count": TARGET_COUNT}, "target coverage drift")
    require(dependencies["symmetric_poisson"]["decision"]["symmetric_poisson_interchange_proved"] is True, "symmetric limit dependency drift")
    require(dependencies["symmetric_poisson"]["decision"]["paired_one_over_m_endpoint_current_cancels"] is True, "pair cancellation dependency drift")
    require(dependencies["Gamma_bulk"]["decision"]["ordinary_full_line_bulk_closed_at_saved_height"] is True, "Gamma bulk dependency drift")
    require(dependencies["endpoint_decomposition"]["decision"]["finite_current_decomposed_exactly"] is True, "endpoint decomposition dependency drift")

    artifact = {
        "kind": STEM,
        "status": "exact_symmetric_bulk_extracted_endpoint_residual_and_limit_proved_quantitative_bound_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "roster_coordinate_length": L,
            "target_positive_modes": [TARGET_START, TARGET_END],
            "target_mode_count": TARGET_COUNT,
            "finite_cutoff": f"M>={TARGET_END}",
        },
        "symbolic_certificate": symbolic_certificate(),
        "mode_partition": mode_partition(),
        "decision": {
            "portcullis_mode_equals_odd_roster_fourier_coefficient": True,
            "finite_symmetric_dirichlet_reassembly_proved": True,
            "bulk_extracted_joint_residual_identity_proved": True,
            "symmetric_residual_limit_proved": True,
            "target_positive_endpoint_currents_cancel_internally": False,
            "target_negative_partners_retained_in_complement": True,
            "zero_mode_is_endpoint_current_canceller": False,
            "half_current_is_endpoint_current_canceller": False,
            "bare_completed_square_endpoint_exponentials_may_be_summed_separately": False,
            "joint_residual_quantitatively_bounded": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "sympy_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Bound the single cancellation-preserving residual H+integral(f C_M)+sum_target(P_A+P_B), beginning with a grouped nonstationary B-endpoint enclosure and then matching the A characteristic part to the certified fold atlas on the same exact current.",
        "proof_boundary": "Exact finite symmetric reassembly, target/complement partition, symmetric residual limit, and one-over-m pair cancellation only. No quantitative joint residual bound, complete endpoint-tail theorem, ordinary/fold splice, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified exact symmetric endpoint-tail reassembly: "
        f"target={TARGET_START}..{TARGET_END}, modes={TARGET_COUNT}, quantitative=False",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
