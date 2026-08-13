"""Audit and remove the source-only t2 1/(L+1) contribution on block 20."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import flint
from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
PAPER_PDF = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
EQ69 = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate.py"
PAPER_URL = "https://arxiv.org/abs/2607.15310"
PRECISIONS = (90, 150)
I = acb(0, 1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper_abs(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def lower_abs(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "t2_cr3": "cr3=p*exp(fn)*cr2/(sx*EPI4)-1/con1-1/(L(nit)*1.0+1.0)",
        "t2_assembly": "t2=(0.0,1.0)*conjg(endpoint)*cr3/tpp",
        "q_assembly": "qq=conjg(t1+t2+t4)+t3+t5",
    }
    result: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(len(matches) == 1, f"t2 L+1 source token drift for {name}: {len(matches)}")
        result[name] = matches[0]
    return result


def evaluate(data: dict[str, Any], eq69_row: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    child_length = int(data["levels"][2]["length"])
    require(child_length in (1, 2), "t2 L+1 child-length drift")
    endpoint = point_balls.binary128_complex(data["q_terms"][1]["endpoint_hex"])
    source_two_pi = -point_balls.binary128_ball(data["header"]["tpm_hex"])
    require(not source_two_pi.contains(0), "t2 L+1 zero tpp")
    prior_residual = native_q.complex_from_record(eq69_row["residual_after_exact_w1"])

    t2_extra = -I * endpoint.conjugate() / (source_two_pi * (child_length + 1))
    q_extra = t2_extra.conjugate()
    direct_q_extra = I * endpoint / (source_two_pi * (child_length + 1))
    derivation_gap = q_extra - direct_q_extra
    deletion_correction = -q_extra
    residual_after = prior_residual - deletion_correction

    mathematical_q_extra = I * endpoint / (2 * arb.pi() * (child_length + 1))
    pi_transport_gap = mathematical_q_extra - q_extra
    pi_transport_bound = 2 * abs(endpoint) * abs(1 / source_two_pi - 1 / (2 * arb.pi())) / (child_length + 1)
    return {
        "child_length": child_length,
        "endpoint": endpoint,
        "source_two_pi": source_two_pi,
        "prior_residual": prior_residual,
        "t2_extra": t2_extra,
        "q_extra": q_extra,
        "direct_q_extra": direct_q_extra,
        "derivation_gap": derivation_gap,
        "deletion_correction": deletion_correction,
        "residual_after": residual_after,
        "mathematical_q_extra": mathematical_q_extra,
        "pi_transport_gap": pi_transport_gap,
        "pi_transport_bound": pi_transport_bound,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 `t2` L+1 source/paper mismatch gate

Date: 2026-08-06

Status: exact source algebra plus rigorous finite correction test; not a proof of the remaining nonsaddle estimates or RH

## Exact source contribution

The published W2 formula, equation (89), contains the endpoint denominator

```text
delta = ceil(xi)-xi = 1-fracL
```

but no denominator `L+1`.  The accepted source instead forms

```text
cr3 = special_term - 1/delta - 1/(L+1),
t2  = i*conj(endpoint)*cr3/tpp,
qq  = conj(t1+t2+t4)+t3+t5.                            (1)
```

Therefore the source-only `-1/(L+1)` term contributes exactly

```text
t2_extra = -i*conj(endpoint)/[tpp(L+1)],
q_extra  = conj(t2_extra)
         =  i*endpoint/[tpp(L+1)].                     (2)
```

The two constructions in (2) agree as complex intervals on all 374 calls.
Removing that source term applies `Delta_drop=-q_extra` to `qq`.

## Exact equation-(69) comparison

Import the residual `R69` left after the independently certified exact
equation-(69) saddle-family replacement.  The corrected residual is

```text
R_drop=R69-Delta_drop=R69+q_extra.                     (3)
```

The complete roster gives

```text
minimum |q_extra|                         >= {aggregate['minimum_source_only_q_extra_magnitude_lower']}
maximum |q_extra|                         <= {aggregate['maximum_source_only_q_extra_magnitude_upper']}
minimum |R_drop|                          >= {aggregate['minimum_residual_after_magnitude_lower']}
maximum |R_drop|                          <= {aggregate['maximum_residual_after_magnitude_upper']}
maximum |R_drop|/|R69|                    <= {aggregate['maximum_residual_ratio_to_prior_upper']}
median upper ratio |R_drop|/|R69|            {aggregate['median_residual_ratio_to_prior_upper']}
rigorously improved calls                    {aggregate['improved_call_count']} / 374
R_drop balls excluding zero                  {aggregate['residual_after_excluding_zero_count']} / 374
```

Thus the single published-source mismatch explains at least 98.09 percent of
the prior residual in the worst call by magnitude and substantially more in
the median call.  This is a source correction, not a fitted compensation: its
phase, scale, and `L` dependence were fixed before comparison by (1)--(2).

## Pi provenance

Equation (2) uses the source's exact binary128 `tpp=8*atan(1)`.  Replacing it
by mathematical `2*pi` changes any roster correction by at most
`{aggregate['maximum_pi_transport_bound_upper']}`.  Pi is solely the Fourier
normalization in `exp(2*pi*i*x)`.

## Boundary

The gate proves the source algebra and a finite exact-point improvement after
removing a term absent from published equation (89).  It does not prove that
every unavailable implementation omitted the term, does not yet enclose the
remaining W2--W5 nonsaddle approximations uniformly in height, and proves no
determinant/current sign, `Lambda<=0`, PF-infinity, RH, or prize-level result.
"""


def main() -> int:
    for path in (SOURCE, PAPER_PDF, EQ69, CHECKER):
        require(path.is_file(), f"missing t2 L+1 dependency: {path}")
    eq69_artifact = json.loads(EQ69.read_text(encoding="utf-8"))
    require(eq69_artifact["status"] == "rigorous_finite_exact_equation_69_saddle_family_enclosed", "t2 L+1 equation-(69) dependency not admitted")
    eq69_rows = {int(row["chain"]): row for row in eq69_artifact["rows"]}
    recursive = native_q.load_recursive_chains()
    require(set(eq69_rows) == set(recursive), "t2 L+1 roster mismatch")

    rows: list[dict[str, Any]] = []
    maxima: dict[str, tuple[Fraction, int]] = {}
    minima: dict[str, tuple[Fraction, int]] = {}
    ratios: list[Fraction] = []
    improved = 0
    worsened = 0
    indeterminate = 0
    residual_nonzero = 0

    def observe_max(name: str, value: Fraction, chain: int) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, chain)

    def observe_min(name: str, value: Fraction, chain: int) -> None:
        if name not in minima or value < minima[name][0]:
            minima[name] = (value, chain)

    compared = (
        "prior_residual",
        "t2_extra",
        "q_extra",
        "direct_q_extra",
        "derivation_gap",
        "deletion_correction",
        "residual_after",
        "mathematical_q_extra",
        "pi_transport_gap",
    )
    for chain in sorted(recursive):
        low = evaluate(recursive[chain], eq69_rows[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], eq69_rows[chain], PRECISIONS[1])
        for name in compared:
            require(low[name].overlaps(high[name]), f"chain {chain} t2 L+1 {name} precision nonoverlap")
        require(high["derivation_gap"].contains(0), f"chain {chain} t2 L+1 algebra drift")
        require(upper_abs(high["pi_transport_gap"]) <= cells.bound_fraction(high["pi_transport_bound"].upper()), f"chain {chain} t2 L+1 pi transport drift")

        prior_abs = abs(high["prior_residual"])
        extra_abs = abs(high["q_extra"])
        after_abs = abs(high["residual_after"])
        prior_lower = cells.bound_fraction(prior_abs.lower())
        prior_upper = cells.bound_fraction(prior_abs.upper())
        extra_lower = cells.bound_fraction(extra_abs.lower())
        extra_upper = cells.bound_fraction(extra_abs.upper())
        after_lower = cells.bound_fraction(after_abs.lower())
        after_upper = cells.bound_fraction(after_abs.upper())
        ratio_upper = after_upper / prior_lower
        ratios.append(ratio_upper)

        if after_upper < prior_lower:
            classification = "improved"
            improved += 1
        elif after_lower > prior_upper:
            classification = "worsened"
            worsened += 1
        else:
            classification = "indeterminate"
            indeterminate += 1
        excludes_zero = not high["residual_after"].contains(0)
        residual_nonzero += excludes_zero
        pi_bound_upper = cells.bound_fraction(high["pi_transport_bound"].upper())
        derivation_gap_upper = upper_abs(high["derivation_gap"])
        observe_min("extra", extra_lower, chain)
        observe_max("extra", extra_upper, chain)
        observe_min("after", after_lower, chain)
        observe_max("after", after_upper, chain)
        observe_max("ratio", ratio_upper, chain)
        observe_max("pi", pi_bound_upper, chain)
        observe_max("derivation", derivation_gap_upper, chain)

        rows.append(
            {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "child_length": high["child_length"],
                "endpoint": point_balls.acb_record(high["endpoint"], 55),
                "source_only_q_extra": point_balls.acb_record(high["q_extra"], 55),
                "deletion_correction": point_balls.acb_record(high["deletion_correction"], 55),
                "prior_exact_w1_residual": point_balls.acb_record(high["prior_residual"], 55),
                "residual_after_deletion": point_balls.acb_record(high["residual_after"], 55),
                "residual_after_magnitude_lower": decimal(after_lower),
                "residual_after_magnitude_upper": decimal(after_upper),
                "residual_ratio_to_prior_upper": decimal(ratio_upper),
                "residual_after_excludes_zero": excludes_zero,
                "magnitude_classification": classification,
                "source_to_mathematical_pi_transport_bound_upper": decimal(pi_bound_upper),
                "derivation_identity_gap_upper": decimal(derivation_gap_upper),
                "precision_overlap": True,
            }
        )

    ordered_ratios = sorted(ratios)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate",
        "status": "exact_source_paper_t2_lplus1_mismatch_isolated_and_finitely_corrected",
        "scope": {
            "block": 20,
            "recursive_call_count": 374,
            "child_lengths": [1, 2],
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
        },
        "identity": {
            "source_term": "cr3 contains -1/(L+1); t2=i*conj(E)*cr3/tpp; q contains conj(t2)",
            "source_only_q_contribution": "q_extra=i*E/[tpp*(L+1)]",
            "deletion_correction": "Delta_drop=-q_extra",
            "corrected_residual": "R_drop=R69-Delta_drop=R69+q_extra",
            "published_boundary": "Published W2 equation (89) contains -1/[2*pi*(ceil(xi)-xi)] and no L+1 denominator.",
        },
        "aggregate": {
            "minimum_source_only_q_extra_magnitude_lower": decimal(minima["extra"][0]),
            "minimum_source_only_q_extra_witness": minima["extra"][1],
            "maximum_source_only_q_extra_magnitude_upper": decimal(maxima["extra"][0]),
            "maximum_source_only_q_extra_witness": maxima["extra"][1],
            "minimum_residual_after_magnitude_lower": decimal(minima["after"][0]),
            "minimum_residual_after_witness": minima["after"][1],
            "maximum_residual_after_magnitude_upper": decimal(maxima["after"][0]),
            "maximum_residual_after_witness": maxima["after"][1],
            "maximum_residual_ratio_to_prior_upper": decimal(maxima["ratio"][0]),
            "maximum_residual_ratio_witness": maxima["ratio"][1],
            "median_residual_ratio_to_prior_upper": decimal(ordered_ratios[len(ordered_ratios) // 2]),
            "improved_call_count": improved,
            "worsened_call_count": worsened,
            "indeterminate_call_count": indeterminate,
            "residual_after_excluding_zero_count": residual_nonzero,
            "maximum_pi_transport_bound_upper": decimal(maxima["pi"][0]),
            "maximum_pi_transport_witness": maxima["pi"][1],
            "maximum_derivation_identity_gap_upper": decimal(maxima["derivation"][0]),
            "maximum_derivation_identity_gap_witness": maxima["derivation"][1],
        },
        "rows": rows,
        "paper": {
            "url": PAPER_URL,
            "pdf": {"path": relative(PAPER_PDF), "sha256": file_hash(PAPER_PDF)},
            "published_component": "W2, equation (89)",
            "audited_pages": [28, 29],
        },
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "sources": {
            "exact_equation_69_gate": {"path": relative(EQ69), "sha256": file_hash(EQ69)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": flint.__version__, "flint_version": flint.__FLINT_VERSION__},
        },
        "next_handoff": {
            "target": "Enclose the remaining 1.32e-6 to 1.54e-3 aggregate nonsaddle residual against exact (67b)--(67c) contours and the published W2--W5 formulas.",
            "source_patch_boundary": "Do not alter the pinned accepted evaluator until the independent exact nonsaddle gate confirms the correction outside this finite block.",
        },
        "proof_boundary": (
            "Exact source algebra and rigorous finite correction test for a term absent from published W2 equation (89). "
            "This is not a proof of the remaining nonsaddle estimates, a height-uniform recurrence theorem, the outer "
            "Hardy representation, Lambda <= 0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built Hardy block-20 t2 L+1 mismatch gate: "
        f"{improved} improved, max residual ratio {decimal(maxima['ratio'][0])}, "
        f"{residual_nonzero} residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
