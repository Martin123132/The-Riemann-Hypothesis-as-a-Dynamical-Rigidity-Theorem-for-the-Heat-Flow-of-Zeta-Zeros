#!/usr/bin/env python3
"""Prove that the commented cubic ecor phase is already present in gc."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_ecor_redundancy_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_ecor_redundancy_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_ecor_redundancy_gate.py"
PAPER_URL = "https://arxiv.org/abs/2607.15310"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def source_locations(lines: list[str]) -> dict[str, int | list[int]]:
    tokens = {
        "cubic_saddle_displacements": ("wm=3*phicoeff(3,nit)*(con1**2)/(xr(nit)**3)", 2),
        "saddle_locations": ("c=con1/xr(nit)-wm", 2),
        "phase_seeds": ("gc=i*c", 2),
        "full_phase_polynomials": ("gc=gc-phicoeff(j,nit)*(c**j)", 2),
        "phase_exponentials": ("cr1=exp(cr1)", 2),
        "commented_ecor_formula": ("ecor=9*p*(phicoeff(3,nit)**2)/(32.0*(phicoeff(2,nit)**5))", 1),
        "zero_ecor_assignment": ("ecor=0.0", 1),
        "ecor_comment": ("Note the exp(I*ecor) correction", 1),
    }
    locations: dict[str, int | list[int]] = {}
    for name, (token, count) in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(len(matches) == count, f"source token count drift for {name}: {len(matches)}")
        locations[name] = matches[0] if count == 1 else matches
    return locations


def executable_ecor_lines(lines: list[str]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for number, raw in enumerate(lines, 1):
        code = raw.split("!", 1)[0]
        if re.search(r"\becor\b", code, flags=re.IGNORECASE):
            records.append({"line": number, "code": code.strip()})
    return records


def symbolic_identity() -> dict[str, str]:
    a, b, x, phi2, pi = sp.symbols("a b x phi2 pi", nonzero=True)
    c_source = a / x - 3 * b * a**2 / x**3
    g_source = sp.factor(sp.expand(a * c_source - x * c_source**2 / 2 - b * c_source**3))
    expected = (
        a**2 / (2 * x)
        - b * a**3 / x**3
        + sp.Rational(9, 2) * b**2 * a**4 / x**5
        - 27 * b**3 * a**5 / x**7
        + 27 * b**4 * a**6 / x**9
    )
    require(sp.simplify(g_source - expected) == 0, "source cubic phase expansion drift")

    phase = sp.expand(2 * pi * expected.subs(x, 2 * phi2))
    quartic = sp.expand(phase).coeff(a, 4).coeff(b, 2)
    commented = 9 * pi / (32 * phi2**5)
    require(sp.simplify(quartic - commented) == 0, "ecor coefficient mismatch")
    doubled = sp.simplify(quartic + commented)
    require(sp.simplify(doubled - 2 * commented) == 0, "double-count identity drift")

    return {
        "source_saddle_location": sp.sstr(c_source),
        "source_gc_expansion": sp.sstr(expected),
        "phase_quartic_coefficient": sp.sstr(quartic),
        "commented_ecor_coefficient": sp.sstr(commented),
        "coefficient_after_hypothetical_extra_ecor": sp.sstr(doubled),
        "quartic_match_identity": sp.sstr(sp.simplify(quartic - commented)),
    }


def build_note(artifact: dict[str, Any]) -> str:
    locations = artifact["source"]["locations"]
    identity = artifact["identity"]
    return f"""# Hardy block-20 `ecor` redundancy gate

Date: 2026-08-06

Status: exact source algebra; the commented leading cubic phase correction is already present in `gc`

## Source calculation

Write

```text
a = i-phi1,  x = 2 phi2,  b = phi3.
```

For a cubic parent the accepted source forms the saddle approximation

```text
c_src = a/x - 3 b a^2/x^3
```

in both endpoint loops at lines {locations['cubic_saddle_displacements']} and
{locations['saddle_locations']}.  It then evaluates the full cubic Legendre
phase `g(c)=a c-x c^2/2-b c^3` at that `c` (source locations
{locations['phase_seeds']}--{locations['phase_exponentials']}).  Exact expansion gives

```text
g(c_src) = a^2/(2x) - b a^3/x^3 + 9 b^2 a^4/(2x^5)
           - 27 b^3 a^5/x^7 + 27 b^4 a^6/x^9.          (1)
```

Since `x=2 phi2` and the source phase is `2*pi*g`, the coefficient of
`b^2 a^4` already present in (1) is

```text
9*pi/(32*phi2^5).                                      (2)
```

This is exactly the coefficient in the commented `ecor` assignment at line
{locations['commented_ecor_formula']}.  Adding a separate factor
`exp(i*ecor*a^4)` to the current `gc` calculation would therefore double the
coefficient in (2), not restore an absent term.

The source assigns `ecor=0` at line {locations['zero_ecor_assignment']}, but
`ecor` has no executable use after that assignment.  Its only executable
occurrences are the declaration and zero assignment.  The nearby statement
about `exp(I*ecor)` is a stale comment relative to this source body.

## Boundary

This exact identity removes the displayed `ecor` coefficient from the list of
plausible missing leading corrections for the accepted evaluator.  It does
not prove that the saddle approximation is exact, control terms beyond the
source `c` model, reconstruct the unavailable Maple implementation, bound the
contour or finite-`ip` remainders, or establish a height-uniform recurrence,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    for path in (SOURCE, CHECKER):
        require(path.is_file(), f"missing ecor dependency: {path}")
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    locations = source_locations(lines)
    executable = executable_ecor_lines(lines)
    require(len(executable) == 2, "unexpected executable ecor use")
    require("real(dp)" in executable[0]["code"].lower(), "ecor declaration drift")
    require(executable[1]["code"].replace(" ", "").lower() == "ecor=0.0", "ecor assignment drift")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_ecor_redundancy_gate",
        "status": "exact_displayed_ecor_coefficient_already_present_in_source_gc_phase",
        "scope": {
            "phase_degree": 3,
            "source_hierarchy_level": 1,
            "physical_block": 20,
            "identity_type": "exact_symbolic_rational_function",
        },
        "identity": symbolic_identity(),
        "source_audit": {
            "executable_ecor_occurrences": executable,
            "executable_ecor_occurrence_count": len(executable),
            "ecor_used_in_phase_expression": False,
            "classification": "stale_comment_or_prior_formula; enabling it in the current full-gc path would double-count the leading b^2*a^4 phase term",
        },
        "source": {
            "path": relative(SOURCE),
            "sha256": file_hash(SOURCE),
            "locations": locations,
        },
        "paper": {
            "url": PAPER_URL,
            "exact_through": "equation (69)",
            "saddle_approximation_begins": "equation (72)",
            "endpoint_component": "W1, equation (81)",
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "sympy_version": sp.__version__,
        },
        "next_handoff": {
            "remove_candidate": "Do not add the commented ecor coefficient to the current cr1=exp(i*tpp*gc) path.",
            "live_target": "Evaluate the actually omitted finite-ip endpoint terms and then bound the difference between the displayed W1 approximation and the exact contour integrals.",
        },
        "proof_boundary": (
            "Exact algebra for the displayed cubic source phase only. It neither proves the saddle or contour approximation "
            "nor reconstructs the unavailable Maple code, controls finite-ip tails, supplies a uniform recurrence bound, "
            "or establishes Lambda <= 0, PF-infinity, RH, or a prize-level theorem."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print("built Hardy block-20 ecor redundancy gate: exact quartic match, 0 executable phase uses")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
