#!/usr/bin/env python3
"""Embed the translated A-face channel in the common post-A remainder."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "finite_regulator_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "extended_notch_normal_form": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate.json",
    "normal_tangential_split": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate.json",
    "midpoint_Peano_partition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate.json",
    "six_current_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.json",
    "complete_A_face_assembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_complete_absolute_assembly_gate.json",
}

T_HEIGHT = 10_000_000_000
A = 159_577
B = 5_122_421
FIRST_TARGET = 622
OLD_LAST = 39_894
FIRST_SHIFTED = 39_895
LAST_SHIFTED = 39_936
LEFT_EDGE = Fraction(79_789, 2)
RIGHT_EDGE = Fraction(79_873, 2)
R_AFTER_A_TARGET = Decimal("0.0368147039947")
A_FACE_THRESHOLD = Decimal("0.00364")
NON_A_SUFFICIENT_TARGET = R_AFTER_A_TARGET - A_FACE_THRESHOLD


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


def exact_ownership_certificate(finite: dict[str, Any]) -> dict[str, Any]:
    target = set(range(FIRST_TARGET, OLD_LAST + 1))
    a_window = set(range(39_853, LAST_SHIFTED + 1))
    extended = set(range(FIRST_TARGET, LAST_SHIFTED + 1))
    shifted = set(range(FIRST_SHIFTED, LAST_SHIFTED + 1))

    require(target | a_window == extended, "extended-notch union drift")
    require(extended - target == shifted, "translated block set difference drift")
    require(len(shifted) == 42, "translated block cardinality drift")
    require(RIGHT_EDGE - LEFT_EDGE == 42, "translated strip length drift")
    require(
        all(Fraction(m, 1) == LEFT_EDGE + Fraction(2 * j + 1, 2) for j, m in enumerate(sorted(shifted))),
        "translated modes are not the strip midpoints",
    )
    require(B > LAST_SHIFTED, "common cutoff does not contain the translated block")

    rows = finite["finite_equivalence_certificate"]["positive_class_rows"]
    shifted_rows = [row for row in rows if row["sector_id"] == "A_window_outer_pairs"]
    require(len(shifted_rows) == 1, "translated coefficient row missing or duplicated")
    row = shifted_rows[0]
    require(row["mode_range"] == [FIRST_SHIFTED, LAST_SHIFTED], "translated row support drift")
    require(row["count"] == 42, "translated coefficient count drift")
    require(row["joined_summand"] == "P_-m+B_m+B_-m", "translated coefficient ownership drift")
    require(
        row["post_A_coefficients"]
        == {"A_minus": 0, "A_plus": 0, "B_minus": 1, "B_plus": 1, "P_minus": 1, "P_plus": 0},
        "translated coefficient vector drift",
    )

    # NT1 and NT2 both have coefficient -1.  The oriented channel is the
    # primal block minus the translated strip, hence (-block)-(-strip).
    block_coefficient = -1
    translated_strip_inside_bracket = -1
    channel_strip_coefficient = -translated_strip_inside_bracket
    require((block_coefficient, channel_strip_coefficient) == (-1, 1), "A-face orientation drift")

    return {
        "old_target": [FIRST_TARGET, OLD_LAST],
        "A_window": [39_853, LAST_SHIFTED],
        "extended_target": [FIRST_TARGET, LAST_SHIFTED],
        "translated_block": [FIRST_SHIFTED, LAST_SHIFTED],
        "translated_block_count": len(shifted),
        "translated_strip": ["79789/2", "79873/2"],
        "translated_strip_length": "42",
        "midpoint_roster_verified": True,
        "common_cutoff_contains_channel": f"M>=B={B}>39936",
        "unchanged_six_class_row": row,
        "oriented_channel_coefficients": {
            "finite_primal_block": block_coefficient,
            "translated_continuous_strip": channel_strip_coefficient,
        },
        "ownership_rule": "The coefficient mask is unchanged. The analytic channel pairs the existing negative finite block with the negative translated strip and charges their signed difference once.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    budget = artifact["budget_certificate"]
    return f"""# Fixed-regulator embedding of the translated A-face channel

Date: 2026-08-24

Status: exact saved-height channel embedding and complement definition;
the non-A complement is not bounded

Let

```text
T={{622,...,39894}},       U={{622,...,39936}},
Q=U minus T={{39895,...,39936}}.                    (EM1)
```

The post-A finite-regulator remainder `mathfrak R_A,(M,epsilon)` is the
common object of Section 11.443.  Sections 11.449--11.453 prove that, for
every `epsilon>0`, translating the upper notch edge from `39894.5` to
`39936.5` produces exactly

```text
C_U,epsilon-C_T,epsilon
 =-sum_(m in Q) exp(-pi*epsilon*m^2)exp(-2*pi*i*m*s),

h_U,epsilon-h_T,epsilon
 =-q_epsilon(y)1_(39894.5<y<39936.5).               (EM2)
```

After the exact A-face lift and removal of its common unit-modulus phase,
the paired edge-translation channel is therefore

```text
E_(42,epsilon)(x)
 =[-sum_(m=39895)^39936 f_(epsilon,x)(m)]
   -[-integral_(39894.5)^39936.5 f_(epsilon,x)(y)dy]. (EM3)
```

Thus the finite block has coefficient `-1` and the translated strip has
coefficient `+1`.  Both pieces in (EM3) are one channel.  Neither may be
estimated or charged as a second copy of the 42-mode row.

For every common cutoff `M>=B=5122421`, let
`mathfrak E_(A,42),epsilon` denote the exact lifted channel (EM3), and define

```text
mathfrak R_nonA,(M,epsilon)
 :=mathfrak R_A,(M,epsilon)-mathfrak E_(A,42),epsilon.

mathfrak R_A,(M,epsilon)
 =mathfrak E_(A,42),epsilon
  +mathfrak R_nonA,(M,epsilon).                     (EM4)
```

The first term is independent of `M`: its finite support ends at `39936`
and its continuous strip is compact.  Equation (EM4) is an analytic-channel
decomposition of the existing post-A remainder.  It does **not** delete a
coefficient row.  In particular, modes `39895..39936` retain the unchanged
six-class summand

```text
P_-m+B_m+B_-m.                                      (EM5)
```

The A notch remains removed, and no `P_m`, `A_m`, or `A_-m` atom is
reintroduced.

The Abel-zero physical split is legitimate.  On the finite strip,
`exp(-pi*epsilon*y^2)` tends uniformly to one.  On `0<x<=x_16`, the exact
six-current formula has a factor `x` in every retained term and its rigorous
remainder has a factor `x^6`; hence `G_A(y,x)=O(x)` uniformly on the strip.
The physical density is consequently

```text
x^(-5/4)(1-x)^(-1/4) O(x)=O(x^(-1/4)),              (EM6)
```

which is integrable at zero.  The interval `x_16<=x<=1/2` is compact.
Dominated convergence, followed by linearity of the already fixed common
limit order, gives

```text
I_(A,42)=lim_(epsilon down 0)
          mathcal P_t[mathfrak E_(A,42),epsilon],

R_nonA=lim_(epsilon down 0)lim_(M to infinity)
        mathcal P_t[mathfrak R_nonA,(M,epsilon)],

R_after_A=I_(A,42)+R_nonA.                          (EM7)
```

Section 11.458 proves `|I_(A,42)|<0.00364`.  Therefore the inherited target
`R_after_A<0.0368147039947` follows from the sufficient, deliberately
strong condition

```text
|R_nonA| < {budget['non_A_sufficient_absolute_target']},

because R_after_A
 <=|I_(A,42)|+|R_nonA|<0.0368147039947.             (EM8)
```

This gate derives (EM8); it does not prove it.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: no new occurrence of `pi` is introduced.  Equation (EM2)
inherits the Gaussian Abel/Fourier normalization, (EM3) inherits the
canonical Fresnel/Kummer phase, and (EM6)--(EM8) introduce no `pi`.

Proof boundary: the translated 42-mode A-face channel is embedded exactly in
the saved-height common-regulator `R_after_A`, its physical Abel-zero split is
justified, and the sufficient non-A target is derived.  The non-A complement,
joined `R_after_A`, `R_Dir`, `Q_K-T`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, and every prize-level conclusion remain unproved.
"""


def main() -> int:
    started = time.time()
    ctx.dps = 120
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "an embedding dependency is not passed")

    finite = dependencies["finite_regulator_equivalence"]
    normal = dependencies["extended_notch_normal_form"]
    split = dependencies["normal_tangential_split"]
    peano = dependencies["midpoint_Peano_partition"]
    six = dependencies["six_current_reduction"]
    assembly = dependencies["complete_A_face_assembly"]

    require(finite["scope"]["height"] == T_HEIGHT, "finite-regulator height drift")
    require(finite["decision"]["finite_R_after_A_defined_on_one_common_regulator"] is True, "common regulator missing")
    require(finite["decision"]["common_kernel_equals_six_class_post_A_mask"] is True, "six-class equivalence missing")
    require(normal["decision"]["post_A_positive_bulk_mask_is_one_contiguous_extended_notch"] is True, "extended notch missing")
    require(split["decision"]["edge_translation_equals_explicit_42_mode_block"] is True, "edge/block identity missing")
    require(peano["decision"]["fixed_regulator_42_mode_and_translated_strip_channel_oriented_exactly"] is True, "oriented channel missing")
    require(peano["decision"]["zero_regulator_compact_strip_limit_proved"] is True, "compact-strip Abel passage missing")
    require(six["decision"]["D6_over_x_regular_at_zero_proved"] is True, "lower-endpoint regularity missing")
    require(six["decision"]["six_current_Fresnel_tail_with_finite_remainder_proved"] is True, "tail remainder theorem missing")
    require(assembly["decision"]["complete_translated_42_mode_A_face_absolute_bound_below_0_point_00364_proved"] is True, "complete A-face bound missing")
    require(CHECKER.is_file(), "independent checker missing")

    ownership = exact_ownership_certificate(finite)
    complete_ball = arb(assembly["assembly_certificate"]["complete_translated_A_face_absolute_bound"])
    require(complete_ball.upper() < arb(str(A_FACE_THRESHOLD)), "A-face threshold drift")
    require(NON_A_SUFFICIENT_TARGET == Decimal("0.0331747039947"), "non-A target arithmetic drift")
    require(NON_A_SUFFICIENT_TARGET > Decimal("0.03317"), "non-A display allowance drift")

    artifact = {
        "kind": STEM,
        "status": "exact_fixed_regulator_A_face_embedding_and_physical_non_A_complement_split_certified_non_A_bound_open",
        "passed": True,
        "scope": {
            "height": T_HEIGHT,
            "A": A,
            "common_cutoff": f"integer M>=B={B}",
            "regulator": "epsilon>0, followed by epsilon down to zero after M to infinity",
            "translated_modes": [FIRST_SHIFTED, LAST_SHIFTED],
            "translated_strip": ["79789/2", "79873/2"],
        },
        "ownership_certificate": ownership,
        "fixed_regulator_embedding": {
            "source_remainder": "mathfrak R_A,(M,epsilon) from the unchanged six-class common-regulator kernel",
            "A_face_channel": "mathfrak E_(A,42),epsilon=Lift_A{[-sum_(m=39895)^39936 f_(epsilon,x)(m)]-[-integral_(39894.5)^39936.5 f_(epsilon,x)(y)dy]}",
            "non_A_complement_definition": "mathfrak R_nonA,(M,epsilon)=mathfrak R_A,(M,epsilon)-mathfrak E_(A,42),epsilon",
            "exact_split": "mathfrak R_A,(M,epsilon)=mathfrak E_(A,42),epsilon+mathfrak R_nonA,(M,epsilon)",
            "cutoff_independence": "mathfrak E_(A,42),epsilon is independent of M for every M>=B because its block ends at 39936 and its strip is compact",
            "mask_guard": "No coefficient is deleted or added; the row 39895..39936 remains P_-m+B_m+B_-m.",
        },
        "physical_limit_certificate": {
            "uniform_Abel_strip_limit": "exp(-pi*epsilon*y^2)->1 uniformly for 39894.5<=y<=39936.5",
            "lower_endpoint_domination": "The certified six-current formula is O(x) and its exact remainder is O(x^6), so G_A(y,x)=O(x) uniformly on the finite strip.",
            "physical_majorant": "x^(-5/4)(1-x)^(-1/4)O(x)=O(x^(-1/4)) at x=0; the endpoint layer is compact",
            "limit_split": "R_after_A=I_(A,42)+R_nonA in the inherited M-to-infinity then epsilon-down-to-zero order",
        },
        "budget_certificate": {
            "R_after_A_target": str(R_AFTER_A_TARGET),
            "complete_A_face_absolute_threshold": str(A_FACE_THRESHOLD),
            "non_A_sufficient_absolute_target": str(NON_A_SUFFICIENT_TARGET),
            "triangle_inequality": "R_after_A<=|I_(A,42)|+|R_nonA|",
            "display_allowance": ">0.03317",
        },
        "decision": {
            "fixed_regulator_A_face_channel_embedded_exactly": True,
            "A_face_channel_independent_of_M_after_common_cutoff": True,
            "common_positive_regulator_preserved": True,
            "block_minus_translated_strip_orientation_verified": True,
            "six_class_mask_unchanged": True,
            "block_and_strip_charged_exactly_once": True,
            "A_notch_atoms_not_reintroduced": True,
            "physical_Abel_zero_limit_split_proved": True,
            "complete_A_face_absolute_bound_applies_to_embedded_channel": True,
            "remaining_non_A_sufficient_target_0_point_0331747039947_derived": True,
            "remaining_non_A_bound_proved": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"elapsed_seconds": round(time.time() - started, 3), "workers": 1},
        "next_obligation": "Decompose and estimate the exact non-A complement without changing the six-class mask. A sufficient absolute target is 0.0331747039947; retain all zero, lower-edge, B-owned, negative-bulk, and remote-positive channels until a cancellation-safe partition is proved.",
        "proof_boundary": "Exact fixed-regulator embedding, coefficient ownership, physical Abel-zero split, and sufficient non-A target arithmetic at t=10^10 only. No bound for the non-A complement, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified exact fixed-regulator translated A-face embedding; non-A bound remains open", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
